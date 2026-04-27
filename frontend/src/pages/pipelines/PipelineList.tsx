/**
 * Pipeline List Page — TanStack Query migration.
 */
import type React from 'react'
import { useMemo, useState } from 'react'
import {
  Card,
  Row,
  Col,
  Button,
  Tag,
  Space,
  Select,
  Modal,
  Form,
  Input,
  InputNumber,
  Switch,
  Tooltip,
  Divider,
  Alert,
  Badge,
  Typography,
} from 'antd'
import {
  ThunderboltOutlined,
  ClockCircleOutlined,
  DatabaseOutlined,
  CloudServerOutlined,
  PlayCircleOutlined,
  InfoCircleOutlined,
  TagsOutlined,
  BulbOutlined,
  RocketOutlined,
  ExperimentOutlined,
} from '@ant-design/icons'

import {
  usePipelineTemplates,
  useExecutePipeline,
  useBatchExecutePipeline,
  useProjectList,
  useSampleList,
} from '@/hooks/queries'
import pipelineService from '@/services/pipeline.service'
import { FilterBar, EnhancedEmptyState, PageSkeleton, FadeIn, StaggeredList } from '@/components/common'
import type { FilterConfig } from '@/components/common'
import { useFilters } from '@/hooks'
import { toast, notify } from '@/utils/notification'
import type { PipelineTemplate, ParamSchema } from '@/types/pipeline'

const { Option } = Select
const { Title, Text } = Typography

export const PipelineList: React.FC = () => {
  const { data: templates, isLoading, error, refetch } = usePipelineTemplates({ is_active: true })

  const { filters, setFilter, resetFilters } = useFilters({
    initialFilters: { search: '', category: 'all' },
  })

  const [executeModalOpen, setExecuteModalOpen] = useState(false)
  const [selectedTemplate, setSelectedTemplate] = useState<PipelineTemplate | null>(null)
  const [batchMode, setBatchMode] = useState(false)
  const [recommendLoading, setRecommendLoading] = useState(false)
  const [selectedProjectInForm, setSelectedProjectInForm] = useState<string>('')
  const [form] = Form.useForm()

  const { data: projectList } = useProjectList({ limit: 200 })
  const items = (projectList as any)?.items ?? (projectList as any)?.data ?? []

  const { data: samplesData } = useSampleList(selectedProjectInForm || undefined)
  const samples = (samplesData as any)?.items ?? (samplesData as any) ?? []

  const executeMutation = useExecutePipeline()
  const batchExecuteMutation = useBatchExecutePipeline()

  const filteredTemplates = useMemo(() => {
    const list = (templates as any)?.items ?? []
    let filtered = list

    if (filters.category !== 'all') {
      filtered = filtered.filter((t: PipelineTemplate) => t.category === filters.category)
    }

    if (filters.search) {
      const search = filters.search.toLowerCase()
      filtered = filtered.filter(
        (t: PipelineTemplate) =>
          t.display_name.toLowerCase().includes(search) ||
          t.description?.toLowerCase().includes(search) ||
          t.tags.some((tag) => tag.toLowerCase().includes(search)),
      )
    }

    return filtered as PipelineTemplate[]
  }, [filters, templates])

  const categories = useMemo(() => {
    const list = (templates as any)?.items ?? []
    return Array.from(new Set(list.map((t: PipelineTemplate) => t.category)))
  }, [templates])

  const filterConfigs: FilterConfig[] = [
    {
      type: 'search',
      key: 'search',
      placeholder: 'Search pipelines by name, description, or tags...',
    },
    {
      type: 'select',
      key: 'category',
      label: 'Category',
      options: [
        { label: 'All Categories', value: 'all' },
        ...categories.map((cat: any) => ({ label: cat, value: cat })),
      ],
    },
  ]

  const handleExecute = (template: PipelineTemplate) => {
    setSelectedTemplate(template)
    setBatchMode(false)
    form.setFieldsValue({
      task_name: `${template.display_name} - ${new Date().toLocaleDateString()}`,
      task_name_prefix: template.display_name,
      parameters: template.default_params,
    })
    setExecuteModalOpen(true)
  }

  const handleProjectChange = (projectId: string) => {
    setSelectedProjectInForm(projectId)
  }

  const handleGetRecommendations = async () => {
    if (!selectedTemplate) {
      return
    }

    const loadingToast = toast.loading('Analyzing historical tasks...')
    setRecommendLoading(true)

    try {
      const projectId = form.getFieldValue('project_id')
      const recommendations = await pipelineService.getParameterRecommendations(selectedTemplate.id, projectId)

      form.setFieldsValue({ parameters: recommendations.recommended_params })

      loadingToast()
      toast.success(
        `Parameters updated! ${recommendations.explanation} (Confidence: ${Math.round(
          recommendations.confidence_score * 100,
        )}%)`,
      )
    } catch (error: any) {
      loadingToast()
      toast.error(`Failed to get recommendations: ${error.message || 'Unknown error'}`)
    } finally {
      setRecommendLoading(false)
    }
  }

  const handleExecuteSubmit = async () => {
    if (!selectedTemplate) {
      return
    }

    try {
      const values = await form.validateFields()

      if (batchMode) {
        if (!values.sample_ids || values.sample_ids.length === 0) {
          toast.warning('Please select at least one sample for batch execution')
          return
        }

        const loadingToast = toast.loading('Creating batch tasks...')

        const result = await batchExecuteMutation.mutateAsync({
          template_id: selectedTemplate.id,
          project_id: values.project_id,
          sample_ids: values.sample_ids,
          task_name_prefix: values.task_name_prefix,
          parameters: values.parameters || selectedTemplate.default_params,
        })

        loadingToast()

        if (result.failed_samples.length > 0) {
          toast.warning(`${result.total_tasks} tasks created, ${result.failed_samples.length} failed`)
        } else {
          notify.success('Batch Execution Started', `${result.total_tasks} tasks created successfully!`)
        }
      } else {
        const loadingToast = toast.loading('Starting pipeline execution...')

        await executeMutation.mutateAsync({
          template_id: selectedTemplate.id,
          task_name: values.task_name,
          project_id: values.project_id,
          sample_ids: values.sample_ids || [],
          parameters: values.parameters || selectedTemplate.default_params,
        })

        loadingToast()
        notify.success('Pipeline Execution Started', 'Your pipeline task has been created and will start shortly.')
      }

      setExecuteModalOpen(false)
      form.resetFields()
    } catch (error: any) {
      if (error.message) {
        notify.error('Execution Failed', error.message)
      }
    }
  }

  const renderParamField = (_key: string, schema: ParamSchema) => {
    switch (schema.type) {
      case 'integer':
      case 'number':
        return (
          <InputNumber
            style={{ width: '100%' }}
            min={schema.min}
            max={schema.max}
            step={schema.step || 1}
            placeholder={`Enter ${schema.label}`}
          />
        )
      case 'boolean':
        return <Switch />
      case 'select':
        return (
          <Select placeholder={`Select ${schema.label}`}>
            {schema.options?.map((opt) =>
              typeof opt === 'string' ? (
                <Option key={opt} value={opt}>
                  {opt}
                </Option>
              ) : (
                <Option key={opt.value} value={opt.value}>
                  {opt.label}
                </Option>
              ),
            )}
          </Select>
        )
      default:
        return <Input placeholder={`Enter ${schema.label}`} />
    }
  }

  if (isLoading && !templates) {
    return <PageSkeleton hasHeader hasFilters rows={6} />
  }

  if (error) {
    return (
      <FadeIn>
        <Alert
          message="Error Loading Pipelines"
          description={(error as Error)?.message || 'Failed to load pipeline templates'}
          type="error"
          showIcon
          action={<Button onClick={() => refetch()}>Retry</Button>}
        />
      </FadeIn>
    )
  }

  return (
    <div>
      <FadeIn direction="up" delay={0} duration={300}>
        <Space direction="vertical" size="middle" style={{ width: '100%', marginBottom: 24 }}>
          <Space align="center">
            <RocketOutlined style={{ fontSize: 32, color: 'var(--color-primary)' }} />
            <div>
              <Title level={2} style={{ margin: 0 }}>
                Pipeline Execution
              </Title>
              <Text type="secondary">Execute NGS analysis pipelines on your data</Text>
            </div>
          </Space>

          <FilterBar
            filters={filterConfigs}
            values={filters}
            onFilterChange={setFilter as (key: string, value: any) => void}
            onReset={resetFilters}
          />

          <Space style={{ width: '100%', justifyContent: 'space-between' }}>
            <Badge count={filteredTemplates.length} showZero>
              <Text strong>Available Pipelines</Text>
            </Badge>
            <Text type="secondary">{(templates as any)?.total || 0} total templates</Text>
          </Space>
        </Space>
      </FadeIn>

      {filteredTemplates.length === 0 && !isLoading ? (
        <FadeIn direction="up" delay={100}>
          <EnhancedEmptyState
            type={filters.search || filters.category !== 'all' ? 'noSearchResults' : 'noData'}
            title={filters.search || filters.category !== 'all' ? 'No matching pipelines' : 'No pipelines available'}
            description={
              filters.search || filters.category !== 'all'
                ? 'Try adjusting your search criteria or filters'
                : 'No pipeline templates are currently available'
            }
            action={
              filters.search || filters.category !== 'all'
                ? { text: 'Clear Filters', onClick: resetFilters }
                : undefined
            }
            size="default"
          />
        </FadeIn>
      ) : (
        <StaggeredList staggerDelay={80} baseDelay={100} direction="up">
          <Row gutter={[16, 16]}>
            {filteredTemplates.map((template) => (
              <Col key={template.id} xs={24} sm={12} lg={8}>
                <Card
                  hoverable
                  title={
                    <Space>
                      <ThunderboltOutlined style={{ color: 'var(--color-primary)' }} />
                      {template.display_name}
                    </Space>
                  }
                  extra={<Tag color="blue">{template.category}</Tag>}
                  actions={[
                    <Button
                      key="execute"
                      type="primary"
                      icon={<PlayCircleOutlined />}
                      onClick={() => handleExecute(template)}
                      size="large"
                    >
                      Execute
                    </Button>,
                  ]}
                  style={{ height: '100%', display: 'flex', flexDirection: 'column' }}
                >
                  <p style={{ minHeight: 60, color: '#666', marginBottom: 16 }}>{template.description}</p>

                  <Divider style={{ margin: '12px 0' }} />

                  <Space direction="vertical" size="small" style={{ width: '100%' }}>
                    <Space>
                      <ClockCircleOutlined style={{ color: 'var(--color-primary)' }} />
                      <Text type="secondary">{template.estimated_time || 'N/A'}</Text>
                    </Space>
                    <Space>
                      <DatabaseOutlined style={{ color: 'var(--color-primary)' }} />
                      <Text type="secondary">{template.min_memory_gb || 'N/A'}</Text>
                    </Space>
                    <Space>
                      <CloudServerOutlined style={{ color: 'var(--color-primary)' }} />
                      <Text type="secondary">{template.min_cpu_cores || 'N/A'} cores</Text>
                    </Space>
                    {template.tags.length > 0 && (
                      <Space wrap>
                        <TagsOutlined style={{ color: 'var(--color-primary)' }} />
                        {template.tags.map((tag) => (
                          <Tag key={tag} color="processing">
                            {tag}
                          </Tag>
                        ))}
                      </Space>
                    )}
                  </Space>
                </Card>
              </Col>
            ))}
          </Row>
        </StaggeredList>
      )}

      <Modal
        title={
          <Space>
            <ThunderboltOutlined style={{ color: 'var(--color-primary)' }} />
            <span style={{ fontWeight: 600 }}>Execute Pipeline: {selectedTemplate?.display_name}</span>
          </Space>
        }
        open={executeModalOpen}
        onCancel={() => {
          setExecuteModalOpen(false)
          form.resetFields()
        }}
        onOk={handleExecuteSubmit}
        width={900}
        okText={
          <Space>
            <RocketOutlined />
            {batchMode ? 'Batch Execute' : 'Execute Now'}
          </Space>
        }
        okButtonProps={{ size: 'large' }}
        cancelButtonProps={{ size: 'large' }}
        confirmLoading={executeMutation.isPending || batchExecuteMutation.isPending}
        destroyOnClose
      >
        <Alert
          message="Pipeline Execution"
          description={selectedTemplate?.description}
          type="info"
          showIcon
          icon={<ExperimentOutlined />}
          style={{ marginBottom: 24 }}
        />

        <Form form={form} layout="vertical" size="large">
          <Card size="small" style={{ marginBottom: 16, background: '#f5f5f5' }}>
            <Space direction="vertical" style={{ width: '100%' }}>
              <Space style={{ width: '100%', justifyContent: 'space-between' }}>
                <Text strong>Execution Mode</Text>
                <Switch
                  checked={batchMode}
                  onChange={setBatchMode}
                  checkedChildren="Batch"
                  unCheckedChildren="Single"
                  size="default"
                />
              </Space>
              <Text type="secondary" style={{ fontSize: 13 }}>
                {batchMode ? (
                  <>
                    <RocketOutlined /> Create one task per sample for parallel processing (recommended for large
                    datasets)
                  </>
                ) : (
                  <>
                    <PlayCircleOutlined /> Create a single task for all selected samples
                  </>
                )}
              </Text>
            </Space>
          </Card>

          {batchMode ? (
            <Form.Item
              name="task_name_prefix"
              label="Task Name Prefix"
              rules={[{ required: true, message: 'Please enter task name prefix' }]}
              tooltip="Each task will be named: [Prefix] - [Sample Name]"
            >
              <Input placeholder="Enter task name prefix" />
            </Form.Item>
          ) : (
            <Form.Item
              name="task_name"
              label="Task Name"
              rules={[{ required: true, message: 'Please enter task name' }]}
            >
              <Input placeholder="Enter task name" />
            </Form.Item>
          )}

          <Form.Item name="project_id" label="Project" rules={[{ required: true, message: 'Please select project' }]}>
            <Select placeholder="Select project" onChange={handleProjectChange} showSearch optionFilterProp="children">
              {items.map((p: any) => (
                <Option key={p.id} value={p.id}>
                  {p.name}
                </Option>
              ))}
            </Select>
          </Form.Item>

          <Form.Item
            name="sample_ids"
            label={batchMode ? 'Samples (Required)' : 'Samples (Optional)'}
            rules={batchMode ? [{ required: true, message: 'Please select at least one sample' }] : []}
            tooltip={
              batchMode
                ? 'One task will be created for each selected sample'
                : 'Leave empty to process all samples in the project'
            }
          >
            <Select
              mode="multiple"
              placeholder={batchMode ? 'Select samples for batch processing' : 'Select samples (leave empty for all)'}
              showSearch
              optionFilterProp="children"
            >
              {samples.map((s: any) => (
                <Option key={s.id} value={s.id}>
                  {s.sample_id}
                </Option>
              ))}
            </Select>
          </Form.Item>

          <Divider orientation="left" style={{ marginTop: 24, marginBottom: 16 }}>
            <Space>
              <Text strong style={{ fontSize: 16 }}>
                Pipeline Parameters
              </Text>
            </Space>
          </Divider>

          <Card size="small" style={{ marginBottom: 16, borderColor: '#faad14', background: '#fffbe6' }}>
            <Space style={{ width: '100%', justifyContent: 'space-between' }}>
              <Space>
                <BulbOutlined style={{ color: '#faad14', fontSize: 18 }} />
                <div>
                  <Text strong style={{ color: '#613400' }}>
                    AI Parameter Recommendations
                  </Text>
                  <br />
                  <Text type="secondary" style={{ fontSize: 12 }}>
                    Get optimized parameters based on your successful historical tasks
                  </Text>
                </div>
              </Space>
              <Button
                type="primary"
                icon={<BulbOutlined />}
                onClick={handleGetRecommendations}
                loading={recommendLoading}
                ghost
                style={{ borderColor: '#faad14', color: '#faad14' }}
              >
                Get Recommendations
              </Button>
            </Space>
          </Card>

          {selectedTemplate &&
            Object.entries(selectedTemplate.param_schema).map(([key, schema]) => (
              <Form.Item
                key={key}
                name={['parameters', key]}
                label={
                  <Space>
                    {schema.label}
                    <Tooltip title={`Default: ${schema.default}`}>
                      <InfoCircleOutlined style={{ color: '#999' }} />
                    </Tooltip>
                  </Space>
                }
                initialValue={schema.default}
                valuePropName={schema.type === 'boolean' ? 'checked' : 'value'}
              >
                {renderParamField(key, schema)}
              </Form.Item>
            ))}
        </Form>
      </Modal>
    </div>
  )
}

export default PipelineList
