/**
 * Pipeline List Page - Browse and execute pipelines
 */
import React, { useEffect, useState } from 'react'
import {
  Card,
  Row,
  Col,
  Button,
  Tag,
  Space,
  Select,
  Input,
  Modal,
  Form,
  InputNumber,
  Switch,
  Tooltip,
  Divider,
} from 'antd'
import {
  ThunderboltOutlined,
  ClockCircleOutlined,
  DatabaseOutlined,
  CpuOutlined,
  PlayCircleOutlined,
  InfoCircleOutlined,
  TagsOutlined,
} from '@ant-design/icons'
import { pipelineService } from '../../services/pipeline.service'
import { useProjectStore } from '../../store/projectStore'
import { useSampleStore } from '../../store/sampleStore'
import type { PipelineTemplate, ParamSchema } from '../../types/pipeline'
import { message } from 'antd'

const { Option } = Select
const { Search } = Input
const { TextArea } = Input

export const PipelineList: React.FC = () => {
  const [templates, setTemplates] = useState<PipelineTemplate[]>([])
  const [filteredTemplates, setFilteredTemplates] = useState<PipelineTemplate[]>([])
  const [loading, setLoading] = useState(false)
  const [selectedCategory, setSelectedCategory] = useState<string>('all')
  const [searchText, setSearchText] = useState('')
  const [executeModalOpen, setExecuteModalOpen] = useState(false)
  const [selectedTemplate, setSelectedTemplate] = useState<PipelineTemplate | null>(null)
  const [form] = Form.useForm()

  const { projects, fetchProjects } = useProjectStore()
  const { samples, fetchSamples } = useSampleStore()

  useEffect(() => {
    loadTemplates()
    fetchProjects()
  }, [])

  const loadTemplates = async () => {
    setLoading(true)
    try {
      const response = await pipelineService.getTemplates({ is_active: true })
      setTemplates(response.items)
      setFilteredTemplates(response.items)
    } catch (error: any) {
      message.error(`Failed to load pipelines: ${error.message}`)
    } finally {
      setLoading(false)
    }
  }

  // Filter templates
  useEffect(() => {
    let filtered = templates

    if (selectedCategory !== 'all') {
      filtered = filtered.filter((t) => t.category === selectedCategory)
    }

    if (searchText) {
      const search = searchText.toLowerCase()
      filtered = filtered.filter(
        (t) =>
          t.display_name.toLowerCase().includes(search) ||
          t.description?.toLowerCase().includes(search) ||
          t.tags.some((tag) => tag.toLowerCase().includes(search))
      )
    }

    setFilteredTemplates(filtered)
  }, [selectedCategory, searchText, templates])

  const categories = Array.from(new Set(templates.map((t) => t.category)))

  const handleExecute = (template: PipelineTemplate) => {
    setSelectedTemplate(template)
    form.setFieldsValue({
      task_name: `${template.display_name} - ${new Date().toLocaleDateString()}`,
      parameters: template.default_params,
    })
    setExecuteModalOpen(true)
  }

  const handleProjectChange = (projectId: string) => {
    fetchSamples({ project_id: projectId })
  }

  const handleExecuteSubmit = async () => {
    if (!selectedTemplate) return

    try {
      const values = await form.validateFields()

      await pipelineService.executePipeline({
        template_id: selectedTemplate.id,
        task_name: values.task_name,
        project_id: values.project_id,
        sample_ids: values.sample_ids || [],
        parameters: values.parameters || selectedTemplate.default_params,
      })

      message.success('Pipeline execution started successfully!')
      setExecuteModalOpen(false)
      form.resetFields()
    } catch (error: any) {
      if (error.message) {
        message.error(`Failed to execute pipeline: ${error.message}`)
      }
    }
  }

  const renderParamField = (key: string, schema: ParamSchema) => {
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
            {schema.options?.map((opt) => {
              if (typeof opt === 'string') {
                return (
                  <Option key={opt} value={opt}>
                    {opt}
                  </Option>
                )
              } else {
                return (
                  <Option key={opt.value} value={opt.value}>
                    {opt.label}
                  </Option>
                )
              }
            })}
          </Select>
        )
      default:
        return <Input placeholder={`Enter ${schema.label}`} />
    }
  }

  return (
    <div>
      {/* Filters */}
      <Card style={{ marginBottom: 16 }}>
        <Space style={{ width: '100%', justifyContent: 'space-between' }}>
          <Space>
            <Search
              placeholder="Search pipelines..."
              style={{ width: 300 }}
              value={searchText}
              onChange={(e) => setSearchText(e.target.value)}
              allowClear
            />
            <Select
              value={selectedCategory}
              onChange={setSelectedCategory}
              style={{ width: 200 }}
            >
              <Option value="all">All Categories</Option>
              {categories.map((cat) => (
                <Option key={cat} value={cat}>
                  {cat}
                </Option>
              ))}
            </Select>
          </Space>
        </Space>
      </Card>

      {/* Pipeline Cards */}
      <Row gutter={[16, 16]}>
        {filteredTemplates.map((template) => (
          <Col key={template.id} xs={24} sm={12} lg={8}>
            <Card
              hoverable
              title={
                <Space>
                  <ThunderboltOutlined style={{ color: '#1890ff' }} />
                  {template.display_name}
                </Space>
              }
              extra={<Tag color="blue">{template.category}</Tag>}
              actions={[
                <Button
                  type="primary"
                  icon={<PlayCircleOutlined />}
                  onClick={() => handleExecute(template)}
                >
                  Execute
                </Button>,
              ]}
            >
              <p style={{ minHeight: 60, color: '#666' }}>{template.description}</p>

              <Divider style={{ margin: '12px 0' }} />

              <Space direction="vertical" size="small" style={{ width: '100%' }}>
                <Space>
                  <ClockCircleOutlined />
                  <span>{template.estimated_time || 'N/A'}</span>
                </Space>
                <Space>
                  <DatabaseOutlined />
                  <span>{template.min_memory_gb || 'N/A'}</span>
                </Space>
                <Space>
                  <CpuOutlined />
                  <span>{template.min_cpu_cores || 'N/A'} cores</span>
                </Space>
                {template.tags.length > 0 && (
                  <Space wrap>
                    <TagsOutlined />
                    {template.tags.map((tag) => (
                      <Tag key={tag} color="default">
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

      {/* Execute Modal */}
      <Modal
        title={
          <Space>
            <ThunderboltOutlined />
            Execute Pipeline: {selectedTemplate?.display_name}
          </Space>
        }
        open={executeModalOpen}
        onCancel={() => setExecuteModalOpen(false)}
        onOk={handleExecuteSubmit}
        width={800}
        okText="Execute"
        confirmLoading={loading}
      >
        <Form form={form} layout="vertical">
          <Form.Item
            name="task_name"
            label="Task Name"
            rules={[{ required: true, message: 'Please enter task name' }]}
          >
            <Input placeholder="Enter task name" />
          </Form.Item>

          <Form.Item
            name="project_id"
            label="Project"
            rules={[{ required: true, message: 'Please select project' }]}
          >
            <Select
              placeholder="Select project"
              onChange={handleProjectChange}
              showSearch
              optionFilterProp="children"
            >
              {projects.map((p) => (
                <Option key={p.id} value={p.id}>
                  {p.name}
                </Option>
              ))}
            </Select>
          </Form.Item>

          <Form.Item name="sample_ids" label="Samples (Optional)">
            <Select
              mode="multiple"
              placeholder="Select samples (leave empty for all)"
              showSearch
              optionFilterProp="children"
            >
              {samples.map((s) => (
                <Option key={s.id} value={s.id}>
                  {s.sample_name}
                </Option>
              ))}
            </Select>
          </Form.Item>

          <Divider>Pipeline Parameters</Divider>

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
