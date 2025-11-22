/**
 * Sample List Page - Sample management with CSV import
 */
import React, { useEffect, useState } from 'react'
import {
  Button,
  Space,
  Select,
  Upload,
  Tag,
} from 'antd'
import {
  PlusOutlined,
  UploadOutlined,
  ExperimentOutlined,
} from '@ant-design/icons'
import type { ColumnsType } from 'antd/es/table'
import { useSampleStore } from '../../store/sampleStore'
import { useProjectStore } from '../../store/projectStore'
import { PageHeader, DataTable } from '../../components/common'
import { toast, notifications } from '../../utils/notification'
import type { Sample } from '../../types/sample'
import dayjs from 'dayjs'

const { Option } = Select

export const SampleList: React.FC = () => {
  const { samples, loading, fetchSamples, importFromCSV, deleteSample } = useSampleStore()
  const { projects, fetchProjects } = useProjectStore()
  const [selectedProject, setSelectedProject] = useState<string>('')

  useEffect(() => {
    fetchProjects()
  }, [])

  useEffect(() => {
    if (selectedProject) {
      fetchSamples({ project_id: selectedProject })
    }
  }, [selectedProject])

  const handleCSVImport = async (file: File) => {
    if (!selectedProject) {
      toast.warning('请先选择一个项目')
      return false
    }

    const loadingToast = toast.loading('导入中...')
    try {
      await importFromCSV(selectedProject, file)
      loadingToast()
      notifications.uploadSuccess()
      // Refresh samples list
      fetchSamples({ project_id: selectedProject })
    } catch (error) {
      loadingToast()
      notifications.uploadError()
    }
    return false // Prevent auto upload
  }

  const columns: ColumnsType<Sample> = [
    {
      title: 'Sample Name',
      dataIndex: 'sample_name',
      key: 'sample_name',
      render: (name) => (
        <Space>
          <ExperimentOutlined style={{ color: 'var(--color-success)' }} />
          <span style={{ fontWeight: 500 }}>{name}</span>
        </Space>
      ),
    },
    {
      title: 'Type',
      dataIndex: 'sample_type',
      key: 'sample_type',
      width: 150,
      render: (type) => type && <Tag color="blue">{type}</Tag>,
    },
    {
      title: 'Files',
      dataIndex: 'file_count',
      key: 'file_count',
      width: 100,
      render: (count) => count || 0,
    },
    {
      title: 'Created',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 150,
      render: (date) => dayjs(date).format('YYYY-MM-DD HH:mm'),
    },
  ]

  return (
    <div>
      <PageHeader
        left={
          <Select
            placeholder="Select Project"
            style={{ width: 300 }}
            value={selectedProject || undefined}
            onChange={setSelectedProject}
          >
            {projects.map((p) => (
              <Option key={p.id} value={p.id}>
                {p.name}
              </Option>
            ))}
          </Select>
        }
        right={
          <>
            <Upload beforeUpload={handleCSVImport} accept=".csv" showUploadList={false}>
              <Button icon={<UploadOutlined />}>Import CSV</Button>
            </Upload>
            <Button type="primary" icon={<PlusOutlined />} disabled={!selectedProject}>
              New Sample
            </Button>
          </>
        }
      />

      <DataTable
        columns={columns}
        dataSource={samples}
        rowKey="id"
        loading={loading}
        pagination={{ pageSize: 20 }}
        emptyText="No Samples"
        emptyDescription="Select a project and add samples or import from CSV"
      />
    </div>
  )
}

export default SampleList
